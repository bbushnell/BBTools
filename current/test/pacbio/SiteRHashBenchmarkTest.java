package test.pacbio;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.HashMap;

import pacbio.SiteR;
import shared.Shared;

/**
 * Finite key controls and a caller-shaped HashMap lookup benchmark for PIC-204.
 * The private complement encoder preserves the superseded mate-1 key policy;
 * exact historical timing requires the old source because production decoding
 * changed with the zero-escape policy.
 * @author Brian Bushnell, Jean
 */
public final class SiteRHashBenchmarkTest{

	/** Runs finite zero-escape controls or one lookup workload.
	 * @param args Optional mode=check|bench, encoder=current|zero|complement,
	 * ids=N, mincpu=N, warmup=N, sparse=t
	 */
	public static void main(String[] args){
		String mode="check", encoderName="current";
		int ids=65536;
		double minCpu=10, warmup=1;
		boolean sparse=false;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected flag=value argument: "+arg);}
			final String key=split[0].toLowerCase();
			final String value=split[1].toLowerCase();
			if("mode".equals(key)){mode=value;}
			else if("encoder".equals(key)){encoderName=value;}
			else if("ids".equals(key)){ids=Integer.parseInt(value);}
			else if("mincpu".equals(key)){minCpu=Double.parseDouble(value);}
			else if("warmup".equals(key)){warmup=Double.parseDouble(value);}
			else if("sparse".equals(key)){sparse=parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}

		if("check".equals(mode)){check();}
		else if("bench".equals(mode)){bench(encoder(encoderName), encoderName, ids, minCpu, warmup, sparse);}
		else{throw new IllegalArgumentException("Unknown mode: "+mode);}
	}

	/** Checks exact zero-escape keys, boundary roundtrips, and mate hash collisions. */
	private static void check(){
		int failures=0;
		failures+=checkKey("zero-id0", 1, 0L, Long.MIN_VALUE, ZERO);
		failures+=checkKey("zero-id1", 1, 1L, -1L, ZERO);
		failures+=checkKey("zero-id2", 1, 2L, -2L, ZERO);
		failures+=checkKey("zero-idmax", 1, Long.MAX_VALUE, Long.MIN_VALUE+1, ZERO);
		failures+=checkKey("complement-id0", 1, 0L, ~0L, COMPLEMENT);
		failures+=checkKey("complement-idmax", 1, Long.MAX_VALUE, ~Long.MAX_VALUE,
				COMPLEMENT);
		for(long id : BOUNDARY_IDS){
			failures+=checkRoundTrip("zero", id, ZERO);
			failures+=checkRoundTrip("complement", id, COMPLEMENT);
			failures+=checkProductionParity(id);
		}
		final int complementCollisions=mateHashCollisions(COMPLEMENT, BOUNDARY_IDS);
		final int zeroCollisions=mateHashCollisions(ZERO, BOUNDARY_IDS);
		failures+=checkScalar("complement_hash_collisions", BOUNDARY_IDS.length,
				complementCollisions);
		failures+=checkScalar("zero_hash_collisions", 0, zeroCollisions);
		if(failures!=0){throw new AssertionError(failures+" SiteR hash control failure(s)");}
	}

	/** Runs one {@code HashMap<Long,SiteR>} construction and repeated paired lookups. */
	private static void bench(Encoder encoder, String encoderName, int ids, double minCpu,
			double warmupSeconds,
			boolean sparse){
		if(ids<1 || ids>Integer.MAX_VALUE/4 || !Double.isFinite(minCpu) ||
				!Double.isFinite(warmupSeconds) || minCpu<=0 || warmupSeconds<0){
			throw new IllegalArgumentException(
					"ids must be in 1..Integer.MAX_VALUE/4, mincpu must be finite/positive, "
					+ "and warmup must be finite/nonnegative");
		}
		final long[] numericIDs=ids(ids, sparse);
		final HashMap<Long, SiteR> map=new HashMap<Long, SiteR>(ids*4);
		for(long id : numericIDs){
			map.put(encoder.encode(id, 0), site(id, 0));
			map.put(encoder.encode(id, 1), site(id, 1));
		}

		final long expected=validateMap(map, numericIDs, encoder);
		final long warmupNanos=secondsToNanos(warmupSeconds);
		final long minCpuNanos=secondsToNanos(minCpu);
		if(warmup(map, numericIDs, encoder, warmupNanos)!=expected){
			throw new AssertionError("Warmup produced a mismatched checksum");
		}
		long checksum=0, passes=0;
		final long wallStart=System.nanoTime(), cpuStart=cpuTime();
		long cpuNanos;
		do{
			final long passChecksum=lookupPass(map, numericIDs, encoder);
			if(passChecksum!=expected){
				throw new AssertionError("Pass checksum "+passChecksum+" != expected "+expected);
			}
			checksum+=passChecksum;
			passes++;
			cpuNanos=cpuTime()-cpuStart;
		}while(cpuNanos<minCpuNanos);
		final long wallNanos=System.nanoTime()-wallStart;
		final long hits=passes*numericIDs.length*2L;
		System.out.println("encoder\tids\tsparse\twarmup_sec\tmincpu_sec\tpasses\t"
				+ "hits\texpected_checksum\tchecksum\twall_nanos\tcpu_nanos");
		System.out.println(encoderName+"\t"+ids+"\t"+sparse+"\t"+warmupSeconds+"\t"+
				minCpu+"\t"+passes+"\t"+hits+"\t"+expected+"\t"+checksum+"\t"+
				wallNanos+"\t"+cpuNanos);
	}

	/** Builds positive dense IDs or a positive odd-multiplier permutation modulo 2^63. */
	private static long[] ids(int ids, boolean sparse){
		final long[] array=new long[ids];
		for(int i=0; i<array.length; i++){
			array[i]=(sparse ? (((i+1)*SPARSE_MULT)&Long.MAX_VALUE) : i+1);
		}
		return array;
	}

	/** Constructs a SiteR payload for one key. */
	private static SiteR site(long id, int pairnum){
		final int start=(int)(id&0x3fffffff);
		return new SiteR(start, start+1, 3, Shared.PLUS, id, pairnum);
	}

	/** Checks one exact encoded key and its logical decoders. */
	private static int checkKey(String name, int pairnum, long id, long expected,
			Encoder encoder){
		final long key=encoder.encode(id, pairnum);
		int failures=0;
		failures+=checkScalar(name+" key", expected, key);
		failures+=checkScalar(name+" id", id, encoder.numericID(key));
		failures+=checkScalar(name+" pair", pairnum, encoder.pairNum(key));
		return failures;
	}

	/** Checks both mates across boundary IDs in the nonnegative numericID domain. */
	private static int checkRoundTrip(String name, long id, Encoder encoder){
		int failures=0;
		for(int pairnum=0; pairnum<2; pairnum++){
			final long key=encoder.encode(id, pairnum);
			failures+=checkScalar(name+" roundtrip id pair="+pairnum+" id="+id,
					id, encoder.numericID(key));
			failures+=checkScalar(name+" roundtrip pair pair="+pairnum+" id="+id,
					pairnum, encoder.pairNum(key));
		}
		return failures;
	}

	/** Checks production keys against the explicit zero-escape reference. */
	private static int checkProductionParity(long id){
		int failures=0;
		for(int pairnum=0; pairnum<2; pairnum++){
			final long expected=ZERO.encode(id, pairnum);
			failures+=checkScalar("production parity pair="+pairnum+" id="+id,
					expected, SiteR.idPairnum(id, pairnum));
		}
		return failures;
	}

	/** Counts mate-0/mate-1 Long.hashCode collisions over the given IDs. */
	private static int mateHashCollisions(Encoder encoder, long[] ids){
		int collisions=0;
		for(long id : ids){
			if(Long.hashCode(encoder.encode(id, 0))==Long.hashCode(encoder.encode(id, 1))){
				collisions++;
			}
		}
		return collisions;
	}

	/** Reports one mismatch and returns one, or returns zero when the values match. */
	private static int checkScalar(String label, long expected, long observed){
		final boolean pass=expected==observed;
		System.out.println(label+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Validates every logical key/payload and returns the expected per-pass checksum. */
	private static long validateMap(HashMap<Long, SiteR> map, long[] ids, Encoder encoder){
		if(map.size()!=ids.length*2){
			throw new AssertionError("Map size "+map.size()+" != expected "+(ids.length*2));
		}
		long checksum=0;
		for(long id : ids){
			checksum+=validateSite(map, encoder, id, 0);
			checksum+=validateSite(map, encoder, id, 1);
		}
		return checksum;
	}

	/** Checks one payload stored at its encoded key and returns its checksum term. */
	private static long validateSite(HashMap<Long, SiteR> map, Encoder encoder, long id,
			int pairnum){
		final SiteR site=map.get(encoder.encode(id, pairnum));
		if(site==null){throw new AssertionError("Missing mate "+pairnum+" id "+id);}
		final int start=(int)(id&0x3fffffff);
		if(site.numericID()!=id || site.pairNum()!=pairnum || site.start!=start){
			throw new AssertionError("Wrong payload for mate "+pairnum+" id "+id+
					": numericID="+site.numericID()+", pairnum="+site.pairNum()+
					", start="+site.start);
		}
		return id+pairnum+start;
	}

	/** Runs warmup passes and returns the final observed logical checksum. */
	private static long warmup(HashMap<Long, SiteR> map, long[] ids, Encoder encoder,
			long warmupNanos){
		final long start=cpuTime();
		long checksum;
		do{
			checksum=lookupPass(map, ids, encoder);
		}while(cpuTime()-start<warmupNanos);
		return checksum;
	}

	/** Performs one mate0/mate1 lookup pass and returns a logical payload checksum. */
	private static long lookupPass(HashMap<Long, SiteR> map, long[] ids, Encoder encoder){
		long checksum=0;
		for(long id : ids){
			SiteR site=map.get(encoder.encode(id, 0));
			if(site==null){throw new AssertionError("Missing mate 0 id "+id);}
			checksum+=site.numericID()+site.pairNum()+site.start;
			site=map.get(encoder.encode(id, 1));
			if(site==null){throw new AssertionError("Missing mate 1 id "+id);}
			checksum+=site.numericID()+site.pairNum()+site.start;
		}
		return checksum;
	}

	/** Returns CPU time for the current benchmark thread. */
	private static long cpuTime(){
		if(!CPU.isCurrentThreadCpuTimeSupported()){
			throw new AssertionError("ThreadMXBean current-thread CPU timing is unsupported");
		}
		if(!CPU.isThreadCpuTimeEnabled()){CPU.setThreadCpuTimeEnabled(true);}
		return CPU.getCurrentThreadCpuTime();
	}

	/** Converts finite seconds to nanoseconds without accepting overflow. */
	private static long secondsToNanos(double seconds){
		assert(Double.isFinite(seconds)) : seconds;
		assert(seconds>=0) : seconds;
		if(seconds==0){return 0;}
		final double nanos=seconds*1000000000.0;
		if(nanos>=Long.MAX_VALUE){
			throw new IllegalArgumentException("Duration too large in seconds: "+seconds);
		}
		return Math.max(1, (long)nanos);
	}

	/** Selects the requested private encoder. */
	private static Encoder encoder(String name){
		if("current".equals(name)){return CURRENT;}
		if("zero".equals(name) || "zero-escape".equals(name)){return ZERO;}
		if("complement".equals(name)){return COMPLEMENT;}
		throw new IllegalArgumentException("Unknown encoder: "+name);
	}

	/** Parses BBTools-style booleans. */
	private static boolean parseBoolean(String value){
		return "t".equals(value) || "true".equals(value) || "1".equals(value);
	}

	/** Small pluggable SiteR key policy. */
	private interface Encoder{
		long encode(long id, int pairnum);
		long numericID(long key);
		int pairNum(long key);
	}

	/** Production zero-escape encoder using -id and Long.MIN_VALUE for mate-1 ID 0. */
	private static final Encoder ZERO=new Encoder(){
		@Override
		public long encode(long id, int pairnum){
			assert(id>=0) : "Site-map keys reserve negative values for mate 1: "+id;
			return (pairnum&1)==0 ? id : (-id)|((id-1)&Long.MIN_VALUE);
		}
		@Override
		public long numericID(long key){return key>=0 ? key : (-key)&Long.MAX_VALUE;}
		@Override
		public int pairNum(long key){return key>=0 ? 0 : 1;}
	};

	/** Actual production SiteR helper. */
	private static final Encoder CURRENT=new Encoder(){
		@Override
		public long encode(long id, int pairnum){return SiteR.idPairnum(id, pairnum);}
		@Override
		public long numericID(long key){return key>=0 ? key : (-key)&Long.MAX_VALUE;}
		@Override
		public int pairNum(long key){return key>=0 ? 0 : 1;}
	};

	/** Superseded one's-complement encoder retained as a historical baseline. */
	private static final Encoder COMPLEMENT=new Encoder(){
		@Override
		public long encode(long id, int pairnum){return (pairnum&1)==0 ? id : ~id;}
		@Override
		public long numericID(long key){return key>=0 ? key : ~key;}
		@Override
		public int pairNum(long key){return key>=0 ? 0 : 1;}
	};

	private static final ThreadMXBean CPU=ManagementFactory.getThreadMXBean();
	private static final long SPARSE_MULT=0x5DEECE66DL;
	private static final long[] BOUNDARY_IDS={0L, 1L, 2L, (1L<<32)-1,
			1L<<32, Long.MAX_VALUE};

	/** Static-only fixture. */
	private SiteRHashBenchmarkTest(){}
}
