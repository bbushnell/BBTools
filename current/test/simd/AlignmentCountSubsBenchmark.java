package test.simd;

import idaligner.IDAlignerStatics;
import shared.Shared;
import simd.Vector;

/** Measures repeated countSubs calls with retained score and offset results.
 * Toggles one exact-match query/reference position per call so repeated results
 * are coupled to changing array contents. SIMD mode searches the 32-offset
 * prealignment window; scalar mode scores the supplied hint as a control.
 * @author Jean */
public final class AlignmentCountSubsBenchmark{

	/** Runs warmup and a fixed amount of timed single-threaded work.
	 * @param args workload=interior|left|right, loops=positive integer, and simd=t|f
	 */
	public static void main(final String[] args){
		String workloadName=WORKLOAD_INTERIOR;
		int loops=65536;
		boolean simd=false;
		for(final String arg : args){
			if(arg.startsWith("workload=")){workloadName=arg.substring(9);}
			else if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("countSubs benchmark requires 256-bit Vector support");
		}

		final Workload workload=workloadFor(workloadName);
		final TestCase[] cases=makeCases(workload);
		final boolean oldSimd=Shared.SIMD;
		try{
			Shared.SIMD=simd;
			final long warmup=repeat(cases, WARMUP_LOOPS);
			final long start=System.nanoTime();
			final long checksum=repeat(cases, loops);
			final long elapsed=System.nanoTime()-start;
			System.out.println("COUNT_SUBS_BENCH workload="+workload.name+" simd="+simd+
					" query_length="+QUERY_LENGTH+" ref_length="+workload.refLength+
					" loops="+loops+" calls="+loops+" wall_ns="+elapsed+
					" checksum="+checksum+" warmup_checksum="+warmup+
					" retained=score+offset");
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Repeats prealignment over a fixed bank of mutable test cases.
	 * @param cases Benchmark cases
	 * @param loops Number of calls
	 * @return Retained checksum of substitution scores and returned offsets
	 */
	private static long repeat(final TestCase[] cases, final int loops){
		final int[] pos=new int[2];
		long checksum=0;
		for(int i=0; i<loops; i++){
			final TestCase c=cases[i&(BANKS-1)];
			final int qIndex=(i+c.expectedOffset)%QUERY_LENGTH;
			final byte oldBase=c.query[qIndex];
			final byte base=(oldBase=='A' ? (byte)'C' : (byte)'A');
			c.query[qIndex]=base;
			c.ref[c.expectedOffset+qIndex]=base;
			pos[0]=c.startOffset;
			pos[1]=0;
			final int score=IDAlignerStatics.countSubs(c.query, c.ref, pos, MAX_SUBS);
			if(score!=0 || pos[0]!=c.expectedOffset){
				throw new AssertionError("countSubs "+c.name+" expected="+c.expectedOffset+
						"/0 actual="+pos[0]+"/"+score);
			}
			checksum=mix(mix(mix(checksum, score), pos[0]), base);
		}
		return checksum;
	}

	/** Builds a fixed bank for the selected workload.
	 * @param workload Workload geometry
	 * @return Benchmark cases
	 */
	private static TestCase[] makeCases(final Workload workload){
		final TestCase[] cases=new TestCase[BANKS];
		for(int bank=0; bank<cases.length; bank++){
			final byte[] query=makeQuery(bank);
			final byte[] ref=makeRef(workload.refLength, bank);
			System.arraycopy(query, 0, ref, workload.expectedOffset, query.length);
			cases[bank]=new TestCase(workload.name+"_"+bank, query, ref,
					workload.startOffset, workload.expectedOffset);
		}
		return cases;
	}

	/** Makes one deterministic query.
	 * @param bank Input bank index
	 * @return Query bases
	 */
	private static byte[] makeQuery(final int bank){
		final byte[] query=new byte[QUERY_LENGTH];
		for(int i=0; i<query.length; i++){query[i]=BASES[(i*7+bank*3)&BASE_MASK];}
		return query;
	}

	/** Makes one deterministic non-matching reference background.
	 * @param length Reference length
	 * @param bank Input bank index
	 * @return Reference bases
	 */
	private static byte[] makeRef(final int length, final int bank){
		final byte[] ref=new byte[length];
		for(int i=0; i<ref.length; i++){ref[i]=BASES[(i*5+bank*11+1)&BASE_MASK];}
		return ref;
	}

	/** Returns the requested benchmark geometry.
	 * @param name Workload name
	 * @return Workload geometry
	 */
	private static Workload workloadFor(final String name){
		if(WORKLOAD_INTERIOR.equals(name)){
			return new Workload(name, 256, 64, 64);
		}
		if(WORKLOAD_LEFT.equals(name)){
			return new Workload(name, 200, 0, 0);
		}
		if(WORKLOAD_RIGHT.equals(name)){
			return new Workload(name, 200, 50, 50);
		}
		throw new IllegalArgumentException("Unknown workload: "+name);
	}

	/** Mixes an int value into the checksum.
	 * @param checksum Previous checksum
	 * @param value Retained value
	 * @return Updated checksum
	 */
	private static long mix(long checksum, final int value){
		checksum^=(value&0xFFFFFFFFL)+0x9E3779B97F4A7C15L;
		checksum*=0xBF58476D1CE4E5B9L;
		return Long.rotateLeft(checksum, 11);
	}

	/** Prevents construction of this standalone benchmark. */
	private AlignmentCountSubsBenchmark(){}

	/** Fixed benchmark geometry. */
	private static final class Workload{
		/** Stores one workload geometry.
		 * @param name_ Workload name
		 * @param refLength_ Reference length
		 * @param startOffset_ Input offset
		 * @param expectedOffset_ Expected returned offset
		 */
		private Workload(final String name_, final int refLength_,
				final int startOffset_, final int expectedOffset_){
			name=name_;
			refLength=refLength_;
			startOffset=startOffset_;
			expectedOffset=expectedOffset_;
		}

		/** Workload name. */
		final String name;
		/** Reference length. */
		final int refLength;
		/** Input offset. */
		final int startOffset;
		/** Expected returned offset. */
		final int expectedOffset;
	}

	/** Fixed mutable benchmark input. */
	private static final class TestCase{
		/** Stores one mutable input pair.
		 * @param name_ Case name
		 * @param query_ Mutable query bases
		 * @param ref_ Mutable reference bases
		 * @param startOffset_ Input offset
		 * @param expectedOffset_ Expected returned offset
		 */
		private TestCase(final String name_, final byte[] query_, final byte[] ref_,
				final int startOffset_, final int expectedOffset_){
			name=name_;
			query=query_;
			ref=ref_;
			startOffset=startOffset_;
			expectedOffset=expectedOffset_;
		}

		/** Case name for diagnostics. */
		final String name;
		/** Mutable query bases. */
		final byte[] query;
		/** Mutable reference bases. */
		final byte[] ref;
		/** Input offset. */
		final int startOffset;
		/** Expected returned offset. */
		final int expectedOffset;
	}

	/** Number of mutable input pairs. */
	private static final int BANKS=128;
	/** Warmup loops run before timed work. */
	private static final int WARMUP_LOOPS=8192;
	/** Query length used by the countSubs boundary repair. */
	private static final int QUERY_LENGTH=150;
	/** Maximum substitutions before countSubs can stop scoring. */
	private static final int MAX_SUBS=150;
	/** Base-table mask. */
	private static final int BASE_MASK=3;
	/** Ordinary read bases. */
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
	/** Workload fully contained in the 32-lane vector window. */
	private static final String WORKLOAD_INTERIOR="interior";
	/** Workload that masks negative reference positions. */
	private static final String WORKLOAD_LEFT="left";
	/** Workload that masks beyond-end reference positions. */
	private static final String WORKLOAD_RIGHT="right";
}
